// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header containing Ipopt solver for general constrained optimization problems
 * \author Simon
 */

#ifndef __IPOPTNONLINEARCONSTRAINTSOLVER_H
#define __IPOPTNONLINEARCONSTRAINTSOLVER_H

#include "ipoptIncludes.h"
#include <goast/Optimization.h>

#ifdef GOAST_WITH_IPOPT


/* use IPOPT to solve a constraint minimization problem of type 
 *  minimize f(x) 
 * over all x s.t. 
 * x_l \leq x \leq x_u
 * g_l(x) \leq g_i(x) \leq g_u(x)
 */
template<typename ConfiguratorType>
class IpoptFirstOrderSolverInterface : public Ipopt::TNLP {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;


  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, VectorType> &_constraintsOp;
  const BaseOp<VectorType, SparseMatrixType> &_constraintsGradientOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;

  unsigned int _numTotalDofs;
  unsigned int _numConstraints;

  const std::vector<RealType> &_x_l, &_x_u;
  const std::vector<RealType> &_g_l, &_g_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptFirstOrderSolverInterface( const BaseOp<VectorType, RealType> &energyOp,
                                  const BaseOp<VectorType, VectorType> &gradientOp,
                                  const BaseOp<VectorType, VectorType> &constraintsOp,
                                  const BaseOp<VectorType, SparseMatrixType> &constraintsGradientOp,
                                  const VectorType &startingPoint,
                                  VectorType &nextIter,
                                  const std::vector<RealType> &x_l = -2.e+19,
                                  const std::vector<RealType> &x_u = 2.e+19,
                                  const std::vector<RealType> &gl = -2.e+19,
                                  const std::vector<RealType> &gu = 2.e+19 )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _constraintsOp( constraintsOp ),
            _constraintsGradientOp( constraintsGradientOp ), _startingPoint( startingPoint ), _nextIter( nextIter ),
            _numTotalDofs( startingPoint.size()), _numConstraints( constraintsOp.getTargetDimension()), _x_l( x_l ),
            _x_u( x_u ), _g_l( gl ), _g_u( gu ), _iter( -1 ) {}


  //returns info about the nlp
  bool get_nlp_info( Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                     IndexStyleEnum &index_style ) override {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _constraintsGradientOp.getNNZ();
    //Nonzeros in hessian: 0 because not implemented
    nnz_h_lag = 0;
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  bool get_bounds_info( Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l,
                        Ipopt::Number *g_u ) override {

    for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l[i];
      x_u[i] = _x_u[i];
    }

    for ( Ipopt::Index i = 0; i < _numConstraints; ++i ) {
      g_l[i] = _g_l[i];
      g_u[i] = _g_u[i];
    }
    return true;
  }

  bool get_starting_point( Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number * /*z_L*/,
                           Ipopt::Number * /*z_U*/, Ipopt::Index /*m*/, bool init_lambda,
                           Ipopt::Number * /*lambda*/ ) override {
    for ( int j = 0; j < _startingPoint.size(); ++j )
      x[j] = _startingPoint[j];
    return true;
  }

  bool eval_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number &obj_value ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    obj_value = _energyOp( v );

    return true;
  }

  bool eval_grad_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number *grad_f ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    VectorType res( v.size());
    _gradientOp.apply( v, res );
    for ( int i = 0; i < n; ++i ) grad_f[i] = res[i];

    return true;
  }

  bool eval_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    VectorType s;
    _constraintsOp.apply( v, s );

    for ( int j = 0; j < m; j++ )
      g[j] = s[j];

    return true;
  }

  bool eval_jac_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Index nele_jac,
                   Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values ) override {

    if ( values == nullptr ) {
      // return the structure of the Jacobian
      // this particular Jacobian is dense
      VectorType v( n );
      for ( int j = 0; j < n; ++j ) v[j] = _startingPoint[j];

      SparseMatrixType res;
      _constraintsGradientOp.apply( v, res );

      Ipopt::Index index = 0;
      for ( int k = 0; k < res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( res, k ); it; ++it ) {
          iRow[index] = it.row();
          jCol[index] = it.col();
          index++;
        }
    }
    else {
      // return the values of the Jacobian of the constraints
      // Convert
      VectorType v( n );
      for ( int i = 0; i < n; ++i ) v[i] = x[i];

      SparseMatrixType res;
      _constraintsGradientOp.apply( v, res );

      Ipopt::Index index = 0;
      for ( int k = 0; k < res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( res, k ); it; ++it ) {
          values[index] = it.value();
          index++;
        }

    }

    return true;
  }

  bool intermediate_callback( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, Ipopt::Number /*obj_value*/,
                              Ipopt::Number /*inf_pr*/, Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,
                              Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/,
                              Ipopt::Number /*d_pr*/, Ipopt::Index /*ls_trials*/, const Ipopt::IpoptData *ip_data,
                              Ipopt::IpoptCalculatedQuantities *ip_cq ) override {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter *tnlpAdapter = nullptr;

    if ( ip_cq != nullptr ) {
      Ipopt::OrigIpoptNLP *origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP *> ( Ipopt::GetRawPtr( ip_cq->GetIpoptNLP()));
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == nullptr ) return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter * > ( GetRawPtr( origNLP->nlp()));
    }
    else
      throw std::invalid_argument( "ip_cq == NULL" );
    auto *x = new double[_startingPoint.size()];   // n == _numDofs
    tnlpAdapter->ResortX( *ip_data->curr()->x(), x );
    delete[] x;

    // Plot solution.
    _iter = iter;

    return true;
  }

  void finalize_solution( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/, const Ipopt::Number *x,
                          const Ipopt::Number * /*z_L*/, const Ipopt::Number * /*z_U*/, Ipopt::Index /*m*/,
                          const Ipopt::Number * /*g*/, const Ipopt::Number * /*lambda*/, Ipopt::Number /*obj_value*/,
                          const Ipopt::IpoptData * /*ip_data*/, Ipopt::IpoptCalculatedQuantities *ip_cq ) override {
    for ( int i = 0; i < _nextIter.size(); ++i ) _nextIter[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }

public :
  RealType getNLPError() const { return _nlpError; }

  int getNumIterations() const { return _iter; }

};


template<typename ConfiguratorType>
class IpoptFirstOrderSolver {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;


  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, VectorType> &_constraintsOp;
  const BaseOp<VectorType, SparseMatrixType> &_constraintsGradientOp;

  const RealType _ipoptTol;
  const int _MaxIterations;
  const std::vector<RealType> &_x_l, &_x_u;
  const std::vector<RealType> &_g_l, &_g_u;
  const int _linearSolverTypeIpopt;
  const int _ipoptPrintLevel;

public:
  IpoptFirstOrderSolver( const BaseOp<VectorType, RealType> &energyOp,
                         const BaseOp<VectorType, VectorType> &gradientOp,
                         const BaseOp<VectorType, VectorType> &constraintsOp,
                         const BaseOp<VectorType, SparseMatrixType> &constraintsGradientOp,
                         const int MaxIterations, const RealType Tolerance,
                         const std::vector<RealType> &x_l,
                         const std::vector<RealType> &x_u,
                         const std::vector<RealType> &gl,
                         const std::vector<RealType> &gu, // \todo fix default parameters
                         const int linearSolverTypeIpopt = 0,
                         const int ipoptPrintLevel = 5 )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _constraintsOp( constraintsOp ),
            _constraintsGradientOp( constraintsGradientOp ),
            _ipoptTol( Tolerance ), _MaxIterations( MaxIterations ), _x_l( x_l ), _x_u( x_u ), _g_l( gl ), _g_u( gu ),
            _linearSolverTypeIpopt( linearSolverTypeIpopt ), _ipoptPrintLevel( ipoptPrintLevel ) {}

  virtual ~IpoptFirstOrderSolver() = default;

  void solve( const VectorType &StartPosition, VectorType &Solution ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptFirstOrderSolverInterface<ConfiguratorType> > tsOpt =
            new IpoptFirstOrderSolverInterface<ConfiguratorType>( _energyOp, _gradientOp,
                                                                  _constraintsOp, _constraintsGradientOp,
                                                                  StartPosition, Solution,
                                                                  _x_l, _x_u, _g_l, _g_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

    // Tolerance
    ipoptApp->Options()->SetNumericValue( "tol", _ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue( "acceptable_tol", _ipoptTol );

    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
    ipoptApp->Options()->SetIntegerValue( "max_iter", _MaxIterations );

    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );

    ipoptApp->RethrowNonIpoptException(true);

//    ipoptApp->Options()->SetStringValue( "mu_strategy", "adaptive" );

    // choose linear solver
    switch ( _linearSolverTypeIpopt ) {
      case 0:
        ipoptApp->Options()->SetStringValue( "linear_solver", "MUMPS" );
        break;
      case 1:
        ipoptApp->Options()->SetStringValue( "linear_solver", "pardiso" );
        break;
      case 27:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma27" );
        break;
      case 57:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma57" );
        break;
      case 77:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma77" );
        break;
      case 86:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma86" );
        break;
      case 97:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma97" );
        break;
      default:
        break;
    }

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;

    // Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    outputIpoptStatus( ipoptStatus, true );
  }

  void testDerivative( const VectorType &StartPosition ) const {
    // Set up masking operators:
    VectorType Solution( StartPosition.size());
    Ipopt::SmartPtr<IpoptFirstOrderSolverInterface<ConfiguratorType> > tsOpt =
            new IpoptFirstOrderSolverInterface<ConfiguratorType>( _energyOp, _gradientOp,
                                                                  _constraintsOp, _constraintsGradientOp,
                                                                  StartPosition, Solution,
                                                                  _x_l, _x_u, _g_l, _g_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

    //Derivative Test
    ipoptApp->Options()->SetStringValue( "derivative_test", "first-order" );
//    ipoptApp->Options()->SetStringValue( "derivative_test_print_all", "yes" );
    ipoptApp->Options()->SetNumericValue( "derivative_test_tol", 1e-4 );
    ipoptApp->Options()->SetNumericValue( "derivative_test_perturbation", 1e-8 );

    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
    ipoptApp->Options()->SetIntegerValue( "max_iter", 1 );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;

    // Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
  }
};

template<typename ConfiguratorType, template <typename C> class CHessianType, bool useTriplets=false>
class IpoptSecondOrderSolverInterface : public Ipopt::TNLP {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;
  typedef std::vector<TripletType> TripletListType;

  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, SparseMatrixType> &_hessianOp;
  const BaseOp<VectorType, VectorType> &_constraintsOp;
  const BaseOp<VectorType, SparseMatrixType> &_constraintsGradientOp;
  const CHessianType<ConfiguratorType> &_constraintsHessianOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;

  unsigned int _numTotalDofs;
  unsigned int _numConstraints;

  const std::vector<RealType> &_x_l, &_x_u;
  const std::vector<RealType> &_g_l, &_g_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptSecondOrderSolverInterface( const BaseOp<VectorType, RealType> &energyOp,
                                   const BaseOp<VectorType, VectorType> &gradientOp,
                                   const BaseOp<VectorType, SparseMatrixType > &hessianOp,
                                   const BaseOp<VectorType, VectorType> &constraintsOp,
                                   const BaseOp<VectorType, SparseMatrixType> &constraintsGradientOp,
                                   const CHessianType<ConfiguratorType> &constraintsHessianOp,
                                   const VectorType &startingPoint,
                                   VectorType &nextIter,
                                   const std::vector<RealType> &x_l = -2.e+19,
                                   const std::vector<RealType> &x_u = 2.e+19,
                                   const std::vector<RealType> &gl = -2.e+19,
                                   const std::vector<RealType> &gu = 2.e+19 )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _hessianOp( hessianOp ),
            _constraintsOp( constraintsOp ), _constraintsGradientOp( constraintsGradientOp ),
            _constraintsHessianOp( constraintsHessianOp ), _startingPoint( startingPoint ), _nextIter( nextIter ),
            _numTotalDofs( startingPoint.size()), _numConstraints( constraintsOp.getTargetDimension()), _x_l( x_l ),
            _x_u( x_u ), _g_l( gl ), _g_u( gu ), _iter( -1 ) {}


  //returns info about the nlp
  bool get_nlp_info( Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                     IndexStyleEnum &index_style ) override {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _constraintsGradientOp.getNNZ();

    //Nonzeros in hessian
    //! \todo find a better way to determine the number of nonzero entries
    nnz_h_lag = 0;

//    VectorType v( n );
//    v.setZero();

    SparseMatrixType energy_res;
    _hessianOp.apply( _startingPoint, energy_res );

    for ( int k = 0; k < energy_res.outerSize(); ++k )
      for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it ) {
        if (it.row() >= it.col()) {
          nnz_h_lag++;
        }

      }

    TensorType constraint_res;
    _constraintsHessianOp.apply( _startingPoint, constraint_res );

    for ( auto &res : constraint_res)
      for ( int k = 0; k < res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( res, k ); it; ++it ) {
          if (it.row() >= it.col()) {
            nnz_h_lag++;
          }
        }
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  bool get_bounds_info( Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l,
                        Ipopt::Number *g_u ) override {

    for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l[i];
      x_u[i] = _x_u[i];
    }

    for ( Ipopt::Index i = 0; i < _numConstraints; ++i ) {
      g_l[i] = _g_l[i];
      g_u[i] = _g_u[i];
    }
    return true;
  }

  bool get_starting_point( Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number * /*z_L*/,
                           Ipopt::Number * /*z_U*/, Ipopt::Index /*m*/, bool init_lambda,
                           Ipopt::Number * /*lambda*/ ) override {
//    init_x == true;
//    init_z == false;
//    init_lambda == false;
    if (_startingPoint.size() != n)
      throw std::length_error("Wrong length of argument");
    for ( int j = 0; j < _startingPoint.size(); ++j )
      x[j] = _startingPoint[j];
    return true;
  }

  bool eval_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number &obj_value ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    obj_value = _energyOp( v );

    return true;
  }

  bool eval_grad_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number *grad_f ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i )  {
      v[i] = x[i];
//      std::cout << "x[" << i << "] : " << x[i] << " -> " << v[i]  << std::endl;
    }

    VectorType res( v.size());
    _gradientOp.apply( v, res );
    for ( int i = 0; i < n; ++i ) {
      grad_f[i] = res[i];
//      std::cout << "grad[" << i << "] : " << res[i] << " -> " << grad_f[i]  << std::endl;
    }


    return true;
  }

  bool eval_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) override {
    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    VectorType s;
    _constraintsOp.apply( v, s );

    for ( int j = 0; j < m; j++ )
      g[j] = s[j];

    return true;
  }

  bool eval_jac_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Index nele_jac,
                   Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values ) override {

    if ( values == nullptr ) {
      // return the structure of the Jacobian
      VectorType v( n );
      v.setZero();

      SparseMatrixType res;
      _constraintsGradientOp.apply( v, res );

      Ipopt::Index index = 0;
      for ( int k = 0; k < res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( res, k ); it; ++it ) {
          iRow[index] = it.row();
          jCol[index] = it.col();
          index++;
        }
    }
    else {
      // return the values of the Jacobian of the constraints
      // Convert
      VectorType v( n );
      for ( int i = 0; i < n; ++i ) v[i] = x[i];

      SparseMatrixType res;
      _constraintsGradientOp.apply( v, res );

      Ipopt::Index index = 0;
      for ( int k = 0; k < res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( res, k ); it; ++it ) {
          values[index] = it.value();
          index++;
        }

    }

    return true;
  }

  bool eval_h( Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                                                            Ipopt::Number obj_factor, Ipopt::Index m,
                                                            const Ipopt::Number *lambda, bool new_lambda,
                                                            Ipopt::Index nele_hess, Ipopt::Index *iRow,
                                                            Ipopt::Index *jCol, Ipopt::Number *values ) override {
    if (useTriplets) {
      if ( values == nullptr ) {
        VectorType v( n );
        v.setZero();
        VectorType l( m );
        l.setZero();

        Ipopt::Index index = 0;

        SparseMatrixType energy_res;
        _hessianOp.apply( v, energy_res );

        for ( int k = 0; k < energy_res.outerSize(); ++k )
          for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it ) {
            if (it.row() >= it.col()) {
//            std::cout << index << "\t" << "-1" << "\t" << it.row() << "\t" << it.col() << std::endl;
              iRow[index] = it.row();
              jCol[index] = it.col();
              index++;
            }

          }

        TripletListType constraint_res;
        _constraintsHessianOp.pushTriplets( v, constraint_res, l );
        //      std::cout << "Structure" << std::endl;
        for ( const auto &trip : constraint_res) {
          if (trip.row() >= trip.col()) {
            //              std::cout << index << "\t" << c << "\t" << it.row() << "\t" << it.col() << std::endl;
            iRow[index] = trip.row();
            jCol[index] = trip.col();
            index++;
          }
        }
      }
      else {
        // Convert
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
        VectorType l ( m ); for ( int i = 0; i < m; ++i ) l[i] = lambda[i];

        Ipopt::Index index = 0;

        SparseMatrixType energy_res;
        _hessianOp.apply( v, energy_res );

        for ( int k = 0; k < energy_res.outerSize(); ++k )
          for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it ) {
            if (it.row() >= it.col()) {
//            std::cout << index << "\t" << "-1" << "\t" << it.row() << "\t" << it.col() << "\t" << it.value() << std::endl;
              values[index] = obj_factor * it.value();
              index++;
            }
          }

        TripletListType constraint_res;
        _constraintsHessianOp.pushTriplets( v, constraint_res, l );
        //      std::cout << "Values" << std::endl;
        for ( const auto &trip : constraint_res) {
          if (trip.row() >= trip.col()) {
            //              std::cout << index << "\t" << c << "\t" << it.row() << "\t" << it.col() << std::endl;
            values[index] = trip.value();
            index++;
          }
        }
      }
    }
    else {
      if ( values == nullptr ) {
        VectorType v( n );
        v.setZero();

        Ipopt::Index index = 0;

        SparseMatrixType energy_res;
        _hessianOp.apply( v, energy_res );

        for ( int k = 0; k < energy_res.outerSize(); ++k )
          for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it ) {
            if ( it.row() >= it.col()) {
//            std::cout << index << "\t" << "-1" << "\t" << it.row() << "\t" << it.col() << std::endl;
              iRow[index] = it.row();
              jCol[index] = it.col();
              index++;
            }

          }

        TensorType constraint_res;
        _constraintsHessianOp.apply( v, constraint_res );

//      std::cout << "Structure" << std::endl;
        for ( int c = 0; c < constraint_res.size(); ++c )
          for ( int k = 0; k < constraint_res[c].outerSize(); ++k )
            for ( typename SparseMatrixType::InnerIterator it( constraint_res[c], k ); it; ++it ) {
              if ( it.row() >= it.col()) {
//              std::cout << index << "\t" << c << "\t" << it.row() << "\t" << it.col() << std::endl;
                iRow[index] = it.row();
                jCol[index] = it.col();
                index++;
              }
            }
      }
      else {
        // Convert
        VectorType v( n );
        for ( int i = 0; i < n; ++i ) v[i] = x[i];

        Ipopt::Index index = 0;

        SparseMatrixType energy_res;
        _hessianOp.apply( v, energy_res );

        for ( int k = 0; k < energy_res.outerSize(); ++k )
          for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it ) {
            if ( it.row() >= it.col()) {
//            std::cout << index << "\t" << "-1" << "\t" << it.row() << "\t" << it.col() << "\t" << it.value() << std::endl;
              values[index] = obj_factor * it.value();
              index++;
            }
          }

        TensorType constraint_res;
        _constraintsHessianOp.apply( v, constraint_res );
//      std::cout << "Values" << std::endl;
        for ( int c = 0; c < constraint_res.size(); ++c )
          for ( int k = 0; k < constraint_res[c].outerSize(); ++k )
            for ( typename SparseMatrixType::InnerIterator it( constraint_res[c], k ); it; ++it ) {
              if ( it.row() >= it.col()) {
//              std::cout << index << "\t" << c << "\t" << it.row() << "\t" << it.col() << "\t" << it.value() << std::endl;
                values[index] = lambda[c] * it.value();
                index++;
              }
            }
      }
    }

    return true;
  }

  bool intermediate_callback( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, Ipopt::Number /*obj_value*/,
                              Ipopt::Number /*inf_pr*/, Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,
                              Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/,
                              Ipopt::Number /*d_pr*/, Ipopt::Index /*ls_trials*/, const Ipopt::IpoptData *ip_data,
                              Ipopt::IpoptCalculatedQuantities *ip_cq ) override {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter *tnlpAdapter = nullptr;

    if ( ip_cq != nullptr ) {
      Ipopt::OrigIpoptNLP *origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP *> ( Ipopt::GetRawPtr( ip_cq->GetIpoptNLP()));
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == nullptr ) return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter * > ( GetRawPtr( origNLP->nlp()));
    }
    else
      throw std::invalid_argument( "ip_cq == NULL" );
    auto *x = new double[_startingPoint.size()];   // n == _numDofs
    tnlpAdapter->ResortX( *ip_data->curr()->x(), x );
    delete[] x;

    // Plot solution.
    _iter = iter;

    return true;
  }

  void finalize_solution( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/, const Ipopt::Number *x,
                          const Ipopt::Number * /*z_L*/, const Ipopt::Number * /*z_U*/, Ipopt::Index /*m*/,
                          const Ipopt::Number * /*g*/, const Ipopt::Number * /*lambda*/, Ipopt::Number /*obj_value*/,
                          const Ipopt::IpoptData * /*ip_data*/, Ipopt::IpoptCalculatedQuantities *ip_cq ) override {
    for ( int i = 0; i < _nextIter.size(); ++i ) _nextIter[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }

public :
  RealType getNLPError() const { return _nlpError; }

  int getNumIterations() const { return _iter; }

};

template<typename ConfiguratorType, template <typename C> class CHessianType, bool useTriplets=false>
class IpoptSecondOrderSolver {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;


  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, SparseMatrixType> &_hessianOp;
  const BaseOp<VectorType, VectorType> &_constraintsOp;
  const BaseOp<VectorType, SparseMatrixType> &_constraintsGradientOp;
  const CHessianType<ConfiguratorType> &_constraintsHessianOp;

  const RealType _ipoptTol;
  const int _MaxIterations;
  const std::vector<RealType> &_x_l, &_x_u;
  const std::vector<RealType> &_g_l, &_g_u;
  const int _linearSolverTypeIpopt;
  const int _ipoptPrintLevel;

public:
  IpoptSecondOrderSolver( const BaseOp<VectorType, RealType> &energyOp,
                          const BaseOp<VectorType, VectorType> &gradientOp,
                          const BaseOp<VectorType, SparseMatrixType> &hessianOp,
                          const BaseOp<VectorType, VectorType> &constraintsOp,
                          const BaseOp<VectorType, SparseMatrixType> &constraintsGradientOp,
                          const CHessianType<ConfiguratorType> &constraintsHessianOp,
                          const int MaxIterations,
                          const RealType Tolerance,
                          const std::vector<RealType> &x_l = -2.e+19, const std::vector<RealType> &x_u = 2.e+19,
                          const std::vector<RealType> &gl = { -2.e+19 },
                          const std::vector<RealType> &gu = { 2.e+19 }, // \todo fix default parameters
                          const int linearSolverTypeIpopt = 0,
                          const int ipoptPrintLevel = 5 )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _hessianOp( hessianOp ), _constraintsOp( constraintsOp ),
            _constraintsGradientOp( constraintsGradientOp ), _constraintsHessianOp( constraintsHessianOp ),
            _ipoptTol( Tolerance ), _MaxIterations( MaxIterations ), _x_l( x_l ), _x_u( x_u ), _g_l( gl ), _g_u( gu ),
            _linearSolverTypeIpopt( linearSolverTypeIpopt ), _ipoptPrintLevel( ipoptPrintLevel ) {}

  virtual ~IpoptSecondOrderSolver() = default;

  void solve( const VectorType &StartPosition, VectorType &Solution ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptSecondOrderSolverInterface<ConfiguratorType, CHessianType, useTriplets> > tsOpt =
            new IpoptSecondOrderSolverInterface<ConfiguratorType, CHessianType, useTriplets>( _energyOp, _gradientOp, _hessianOp,
                                                                                 _constraintsOp, _constraintsGradientOp,
                                                                                 _constraintsHessianOp,
                                                                                 StartPosition, Solution,
                                                                                 _x_l, _x_u, _g_l, _g_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

    // Tolerance
    ipoptApp->Options()->SetNumericValue( "tol", _ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue( "acceptable_tol", _ipoptTol );

    ipoptApp->Options()->SetIntegerValue( "max_iter", _MaxIterations );

    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );

    ipoptApp->RethrowNonIpoptException(true);

//    ipoptApp->Options()->SetNumericValue( "bound_relax_factor", 1e-10 );


//    ipoptApp->Options()->SetStringValue( "mu_strategy", "adaptive" );

    // choose linear solver
    switch ( _linearSolverTypeIpopt ) {
      case 0:
        ipoptApp->Options()->SetStringValue( "linear_solver", "MUMPS" );
        break;
      case 1:
        ipoptApp->Options()->SetStringValue( "linear_solver", "pardiso" );
        break;
      case 2:
        ipoptApp->Options()->SetStringValue( "linear_solver", "wsmp" );
        ipoptApp->Options()->SetIntegerValue( "wsmp_num_threads", 4 );
        break;
      case 27:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma27" );
        break;
      case 57:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma57" );
        break;
      case 77:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma77" );
        break;
      case 86:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma86" );
        break;
      case 97:
        ipoptApp->Options()->SetStringValue( "linear_solver", "ma97" );
        break;
      default:
        break;
    }

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;

    // Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    outputIpoptStatus( ipoptStatus, true );
  }

  void testDerivative( const VectorType &StartPosition ) const {
    // Set up masking operators:
    VectorType Solution( StartPosition );
    Ipopt::SmartPtr<IpoptSecondOrderSolverInterface<ConfiguratorType, CHessianType> > tsOpt =
            new IpoptSecondOrderSolverInterface<ConfiguratorType, CHessianType>( _energyOp, _gradientOp, _hessianOp,
                                                                                 _constraintsOp, _constraintsGradientOp,
                                                                                 _constraintsHessianOp,
                                                                                 StartPosition, Solution,
                                                                                 _x_l, _x_u, _g_l, _g_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

    //Derivative Test
    ipoptApp->Options()->SetStringValue( "derivative_test", "only-second-order" );
    ipoptApp->Options()->SetStringValue( "derivative_test_print_all", "no" );
    ipoptApp->Options()->SetNumericValue( "derivative_test_tol", 1e-3 );
    ipoptApp->Options()->SetNumericValue( "derivative_test_perturbation", 1e-12 );
    ipoptApp->Options()->SetNumericValue( "point_perturbation_radius", 0. );

    ipoptApp->Options()->SetIntegerValue( "max_iter", 0 );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;

    // Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
  }
};

#endif // GOAST_WITH_IPOPT

#endif //__IPOPTNONLINEARCONSTRAINTSOLVER_H