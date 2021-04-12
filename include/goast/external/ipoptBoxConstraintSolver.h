// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __IPOPTBOXCONSTRAINTSOLVER_H
#define __IPOPTBOXCONSTRAINTSOLVER_H

#include "ipoptIncludes.h"

#include <goast/Optimization/optInterface.h>

#ifdef GOAST_WITH_IPOPT

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief IPOPT solver interface to solve a first-order constraint minimization problem with box constraints
//! \author Simon
//! Problem: Minimize f(x) over all \f$ x\in \R^n \f$ s.t. \f$ x_l \leq x_i \leq x_u \f$ for \f$ i =1, \ldots, n \f$.
//! However, we only have the gradient Df of f at hand.
template<typename ConfiguratorType>
class IpoptBoxConstraintFirstOrderSolverInterface : public Ipopt::TNLP {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;


  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;

  const VectorType &_startingPoint;
  VectorType& _solution;

  unsigned int _numTotalDofs;
  unsigned int _numConstraints;

  const std::vector<RealType> &_x_l, &_x_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptBoxConstraintFirstOrderSolverInterface( const BaseOp<VectorType, RealType> &energyOp,
                                               const BaseOp<VectorType, VectorType> &gradientOp,
                                               const VectorType &startingPoint,
                                               VectorType &solution,
                                               const std::vector<RealType> &x_l, const std::vector<RealType> &x_u )
          : _energyOp( energyOp ), _gradientOp( gradientOp ),
            _startingPoint( startingPoint ), _solution( solution ),
            _numTotalDofs( startingPoint.size()), _numConstraints( 0 ),
            _x_l( x_l ), _x_u( x_u ), _iter( -1 ) {}

  //returns info about the nlp
  bool get_nlp_info( Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                     IndexStyleEnum &index_style ) override {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _numConstraints * _numTotalDofs;
    //Nonzeros in hessian: 0 because not implemented
    nnz_h_lag = 0;
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  bool get_bounds_info( Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                        Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u ) override {

    for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l[i];
      x_u[i] = _x_u[i];
    }

    return true;
  }

  bool get_starting_point( Ipopt::Index n, bool init_x, Ipopt::Number *x,
                           bool init_z, Ipopt::Number * /*z_L*/, Ipopt::Number * /*z_U*/,
                           Ipopt::Index /*m*/, bool init_lambda, Ipopt::Number * /*lambda*/ ) override {
    init_x == true;
    init_z == false;
    init_lambda == false;
    for ( int j = 0; j < _numTotalDofs; ++j )
        x[j] = _startingPoint[j];
    return true;
  }

  bool eval_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number &obj_value ) override {
    if( n != _numTotalDofs )
        throw BasicException("IpoptBoxConstraintFirstOrderSolverInterface::eval_f(): argument has wrong size!");

    // Convert
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    obj_value = _energyOp( v );
    return true;

  }

  bool eval_grad_f( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Number *grad_f ) override {
    if( n != _numTotalDofs )
        throw BasicException("IpoptBoxConstraintFirstOrderSolverInterface::eval_grad_f(): argument has wrong size!");

    // Convert
    //! \todo Replace copying by Eigen Map
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    VectorType res( v.size());
    _gradientOp.apply( v, res );
    for ( int i = 0; i < n; ++i )
      grad_f[i] = res[i];

    return true;
  }

  bool eval_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) override {
    return true;
  }

  bool eval_jac_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m,
                   Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values ) override {
    return true;
  }

  bool intermediate_callback( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter,
                              Ipopt::Number /*obj_value*/, Ipopt::Number /*inf_pr*/, Ipopt::Number /*inf_du*/,
                              Ipopt::Number /*mu*/, Ipopt::Number /*d_norm*/,
                              Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/,
                              Ipopt::Number /*d_pr*/, Ipopt::Index /*ls_trials*/,
                              const Ipopt::IpoptData *ip_data,
                              Ipopt::IpoptCalculatedQuantities *ip_cq ) override {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter *tnlpAdapter = nullptr;

    if ( ip_cq != nullptr ) {
      Ipopt::OrigIpoptNLP *origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP *> ( Ipopt::GetRawPtr( ip_cq->GetIpoptNLP()));
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == nullptr )
        return true;
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
                          const Ipopt::Number * /*z_L*/, const Ipopt::Number * /*z_U*/,
                          Ipopt::Index /*m*/, const Ipopt::Number * /*g*/, const Ipopt::Number * /*lambda*/,
                          Ipopt::Number /*obj_value*/, const Ipopt::IpoptData * /*ip_data*/,
                          Ipopt::IpoptCalculatedQuantities *ip_cq ) override {
    if( _solution.size() != _numTotalDofs )
      _solution.resize( _numTotalDofs );
    for ( int i = 0; i < _numTotalDofs; ++i )
      _solution[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }

  RealType getNLPError() const { return _nlpError; }

  int getNumIterations() const { return _iter; }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Solver class for first-order constrained problem with box constraints.
//! \author Simon
//! See documentation of class IpoptBoxConstraintFirstOrderSolverInterface<> above.
template<typename ConfiguratorType>
class IpoptBoxConstraintFirstOrderSolver {
protected:

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;


  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;

  RealType _ipoptTol;
  RealType _MaxIterations;
  const std::vector<RealType> &_x_l, &_x_u;
  const int _linearSolverTypeIpopt;
  const int _ipoptPrintLevel;

public:
  mutable SolverStatus<ConfiguratorType> status;
  IpoptBoxConstraintFirstOrderSolver( const BaseOp<VectorType, RealType> &energyOp,
                                      const BaseOp<VectorType, VectorType> &gradientOp,
                                      const int MaxIterations, const RealType Tolerance,
                                      const std::vector<RealType> &x_l, const std::vector<RealType> &x_u,
                                      const int linearSolverTypeIpopt = 0,
                                      const int ipoptPrintLevel = 5
  )
          : _energyOp( energyOp ), _gradientOp( gradientOp ),
            _ipoptTol( Tolerance ), _MaxIterations( MaxIterations ),
            _x_l( x_l ), _x_u( x_u ),
            _ipoptPrintLevel( ipoptPrintLevel ),
            _linearSolverTypeIpopt( linearSolverTypeIpopt ) {}

  ~IpoptBoxConstraintFirstOrderSolver() = default;

  void solve( const VectorType &StartPosition, VectorType &Solution ) const {

    Solution.resize( StartPosition.size() );

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintFirstOrderSolverInterface<ConfiguratorType>> tsOpt
            = new IpoptBoxConstraintFirstOrderSolverInterface<ConfiguratorType>( _energyOp, _gradientOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

    //Derivative Test
#ifdef IPOPTDERIVATIVETEST
    ipoptApp->Options ()->SetStringValue  ( "derivative_test", "first-order" );
//     ipoptApp->Options ()->SetStringValue  ( "derivative_test_print_all", "yes" );
    ipoptApp->Options ()->SetStringValue  ( "derivative_test_print_all", "no" );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_tol", 1.e-4 );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_perturbation", 1.e-6 );
#endif //IPOPTDERIVATIVETEST

//     ipoptApp->Options()->SetStringValue( "output_file", "ipoptOutputTest.txt" );

    // Tolerance
    ipoptApp->Options()->SetNumericValue( "tol", _ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue( "acceptable_tol", _ipoptTol );

    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
    ipoptApp->Options()->SetIntegerValue( "max_iter", _MaxIterations );
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

    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;

    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );

    outputIpoptStatus( ipoptStatus, true );
  }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! \brief IPOPT solver interface to solve a second-order constrained minimization problem with box constraints
//! \author Simon
//! Problem: Minimize f(x) over all \f$ x\in \R^n \f$ s.t. \f$ x_l \leq x_i \leq x_u \f$ for \f$ i =1, \ldots, n \f$.
//! Here, we have the Hessian D^2f of f at hand.
template<typename ConfiguratorType, bool useTriplets=false>
class IpoptBoxConstraintSecondOrderSolverInterface : public Ipopt::TNLP {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, SparseMatrixType> &_hessianOp;

  const VectorType &_startingPoint;
  VectorType& _solution;

  unsigned int _numTotalDofs;
  unsigned int _numConstraints;

  const std::vector<RealType> &_x_l, &_x_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptBoxConstraintSecondOrderSolverInterface( const BaseOp<VectorType, RealType> &energyOp,
                                                const BaseOp<VectorType, VectorType> &gradientOp,
                                                const BaseOp<VectorType, SparseMatrixType> &hessianOp,
                                                const VectorType &startingPoint, VectorType& Solution,
                                                const std::vector<RealType> &x_l, const std::vector<RealType> &x_u )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _hessianOp(hessianOp),
            _startingPoint( startingPoint ), _solution(Solution),
            _numTotalDofs( startingPoint.size()), _numConstraints( 0 ),
            _x_l( x_l ), _x_u( x_u ), _iter( -1 ) {}

  //returns info about the nlp
  bool get_nlp_info( Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                     IndexStyleEnum &index_style ) override {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _numConstraints * _numTotalDofs;
    //Nonzeros in hessian: 0 because not implemented
    //! \todo find a better way to determine the number of nonzero entries
    nnz_h_lag = 0;

    VectorType v( n );
    v.setZero();

    if (useTriplets) {
      TripletListType energy_res;
      _hessianOp.pushTriplets( _startingPoint, energy_res );

      for ( const auto &trip : energy_res)
        if (trip.row() >= trip.col())
          nnz_h_lag++;

//      std::cout << nnz_h_lag << " nonzeros in hessian" << std::endl;
    }
    else {
      SparseMatrixType energy_res;
      _hessianOp.apply( _startingPoint, energy_res );

      for ( int k = 0; k < energy_res.outerSize(); ++k )
        for ( typename SparseMatrixType::InnerIterator it( energy_res, k ); it; ++it )
          if ( it.row() >= it.col())
            nnz_h_lag++;
    }

    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  bool get_bounds_info( Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                        Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u ) override {

    for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l[i];
      x_u[i] = _x_u[i];
    }

    return true;
  }

  bool get_starting_point( Ipopt::Index n, bool init_x, Ipopt::Number *x,
                           bool init_z, Ipopt::Number * /*z_L*/, Ipopt::Number * /*z_U*/,
                           Ipopt::Index /*m*/, bool init_lambda, Ipopt::Number * /*lambda*/ ) override {
//    init_x == true;
//    init_z == false;
//    init_lambda == false;
    for ( int j = 0; j < _startingPoint.size(); ++j ) x[j] = _startingPoint[j];
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
    //! \todo Replace copying by Eigen Map
    VectorType v( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    VectorType res( v.size());
    _gradientOp.apply( v, res );
    for ( int i = 0; i < n; ++i )
      grad_f[i] = res[i];

    return true;
  }

  bool eval_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) override {
    return true;
  }

  bool eval_jac_g( Ipopt::Index n, const Ipopt::Number *x, bool /*new_x*/, Ipopt::Index m,
                   Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values ) override {
    return true;
  }

  bool intermediate_callback( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter,
                              Ipopt::Number /*obj_value*/, Ipopt::Number /*inf_pr*/, Ipopt::Number /*inf_du*/,
                              Ipopt::Number /*mu*/, Ipopt::Number /*d_norm*/,
                              Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/,
                              Ipopt::Number /*d_pr*/, Ipopt::Index /*ls_trials*/,
                              const Ipopt::IpoptData *ip_data,
                              Ipopt::IpoptCalculatedQuantities *ip_cq ) override {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter *tnlpAdapter = nullptr;

    if ( ip_cq != nullptr ) {
      Ipopt::OrigIpoptNLP *origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP *> ( Ipopt::GetRawPtr( ip_cq->GetIpoptNLP()));
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == nullptr )
        return true;
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

  bool eval_h( Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor,
               Ipopt::Index m, const Ipopt::Number *lambda, bool new_lambda,
               Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values ) override {
    if (useTriplets) {
      if ( values == nullptr ) {
//        VectorType v( n );
//        v.setZero();

        Ipopt::Index index = 0;

        TripletListType energy_res;
        _hessianOp.pushTriplets( _startingPoint, energy_res );

        for ( const auto &trip : energy_res) {
          if (trip.row() >= trip.col()) {
            iRow[index] = trip.row();
            jCol[index] = trip.col();
            index++;
          }
        }
      }
      else {
        // Convert
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
        Ipopt::Index index = 0;

        TripletListType energy_res;
        _hessianOp.pushTriplets( v, energy_res );

        for ( const auto &trip : energy_res) {
          if (trip.row() >= trip.col()) {
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
              values[index] = it.value();
              index++;
            }
          }

      }
    }
    return true;
  }

  void finalize_solution( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/, const Ipopt::Number *x,
                          const Ipopt::Number * /*z_L*/, const Ipopt::Number * /*z_U*/,
                          Ipopt::Index /*m*/, const Ipopt::Number * /*g*/, const Ipopt::Number * /*lambda*/,
                          Ipopt::Number /*obj_value*/, const Ipopt::IpoptData * /*ip_data*/,
                          Ipopt::IpoptCalculatedQuantities *ip_cq ) override {
   if( _solution.size() != _numTotalDofs )
      _solution.resize( _numTotalDofs );
   for ( int i = 0; i < _numTotalDofs; ++i )
      _solution[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }

  RealType getNLPError() const { return _nlpError; }

  int getNumIterations() const { return _iter; }

  const VectorType& getSolution() const {
      return _solution;
  }
};

template<typename ConfiguratorType, bool useTriplets=false>
class IpoptBoxConstraintSecondOrderSolver {
protected:

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

  const BaseOp<VectorType, RealType> &_energyOp;
  const BaseOp<VectorType, VectorType> &_gradientOp;
  const BaseOp<VectorType, SparseMatrixType> &_hessianOp;

  RealType _ipoptTol;
  RealType _MaxIterations;
  const std::vector<RealType> &_x_l, &_x_u;
  const int _linearSolverTypeIpopt;
  const int _ipoptPrintLevel;

public:
  mutable SolverStatus<ConfiguratorType> status;

  IpoptBoxConstraintSecondOrderSolver( const BaseOp<VectorType, RealType> &energyOp,
                                       const BaseOp<VectorType, VectorType> &gradientOp,
                                       const BaseOp<VectorType, SparseMatrixType> &hessianOp,
                                       const int MaxIterations, const RealType Tolerance,
                                       const std::vector<RealType> &x_l, const std::vector<RealType> &x_u,
                                       const int linearSolverTypeIpopt = 0,
                                       const int ipoptPrintLevel = 5
  )
          : _energyOp( energyOp ), _gradientOp( gradientOp ), _hessianOp( hessianOp ),
            _ipoptTol( Tolerance ), _MaxIterations( MaxIterations ),
            _x_l( x_l ), _x_u( x_u ),
            _ipoptPrintLevel( ipoptPrintLevel ),
            _linearSolverTypeIpopt( linearSolverTypeIpopt ) {

//    if (useTriplets && ipoptPrintLevel >=5) {
//      std::cout << "NOTE: IpoptBoxConstraintSecondOrderSolver is using triplets!" << std::endl;
//    }
  }

  ~IpoptBoxConstraintSecondOrderSolver() = default;

  void solve( const VectorType &StartPosition, VectorType &Solution ) const {

    Solution.resize( StartPosition.size() );

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintSecondOrderSolverInterface<ConfiguratorType, useTriplets>> tsOpt
            = new IpoptBoxConstraintSecondOrderSolverInterface<ConfiguratorType, useTriplets>( _energyOp, _gradientOp, _hessianOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = new Ipopt::IpoptApplication();

//     ipoptApp->Options()->SetStringValue( "output_file", "ipoptOutputTest.txt" );

    // Tolerance
    ipoptApp->Options()->SetNumericValue( "tol", _ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue( "acceptable_tol", _ipoptTol );

    ipoptApp->Options()->SetNumericValue( "tiny_step_tol", 1e-10 );


    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetIntegerValue( "max_iter", _MaxIterations );
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

    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    // Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );

    // print output message to console
    outputIpoptStatus( ipoptStatus, true );
  }


};


#endif //GOAST_WITH_IPOPT

#endif //__IPOPTBOXCONSTRAINTSOLVER_H