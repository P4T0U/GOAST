// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for the Byrd-Omojokun SQP method for equality-constrained problems
 * \author Sassen
 */

#ifndef OPTIMIZATION_BYRDOMOJOKUNSQP_H
#define OPTIMIZATION_BYRDOMOJOKUNSQP_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include "optInterface.h"

#include "modifiedDogleg.h"
#include "projectedCG.h"

template<typename ConfiguratorType, template<typename M> class LinearSolver=Eigen::UmfPackLU>
class ByrdOmojokunSQP : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  typedef LinearSolver<MatrixType> LinearSolverType;

  // Epsilon for safe-guarding computation of reduction against numerically instabilities
  const RealType eps_red = 100 * std::numeric_limits<RealType>::epsilon();

  // Functionals
  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const BaseOp<VectorType, MatrixType> &_D2F;

  const BaseOp<VectorType, VectorType> &_G;
  const BaseOp<VectorType, MatrixType> &_DG;
  const BaseOp<VectorType, TensorType> &_D2G;

  // Dimensionality of problem
  const int _n;

  // Linear solver for column-/row-/null-space projection
  mutable LinearSolverType _linearSolver;

  // General parameters
  int _maxIterations;
  bool _quiet;
  const std::vector<int> *_fixedVariables;

  // SQP parameters
  RealType _eps, _tau;
  RealType _initRadius = 10;
  RealType _maxRadius = std::numeric_limits<RealType>::infinity();
  RealType _minRadius = 1e-12;
  RealType _minStepsize = 1e-12;
  RealType _eta = 0.25;
  RealType _penaltyFactor = 0.3; // rho from formula (18.35) on p. 542 in Nocedal & Wright
  RealType _socThreshold = 0.1;

  std::string _meritFunction = "L2";

  int _diagonalPreconditioning = -1;

  // Tangential solver parameters
  std::string _tangentialMethod = "ProjectedCG";
  std::map<std::string, int> _tangentialIntParameters;
  std::map<std::string, RealType> _tangentialRealParameters;
  std::map<std::string, std::string> _tangentialStringParameters;

  // Normal solver parameters
  std::string _normalMethod = "ModifiedDogleg";
  std::map<std::string, int> _normalIntParameters;
  std::map<std::string, RealType> _normalRealParameters;
  std::map<std::string, std::string> _normalStringParameters;

public:
  ByrdOmojokunSQP( const BaseOp<VectorType, RealType> &F, const BaseOp<VectorType, VectorType> &DF,
                   const BaseOp<VectorType, MatrixType> &D2F, const BaseOp<VectorType, VectorType> &G,
                   const BaseOp<VectorType, MatrixType> &DG, const BaseOp<VectorType, TensorType> &D2G,
                   RealType tolerance = 1e-8, int maxIterations = 1000,
                   bool quiet = true ) : _F( F ), _DF( DF ), _D2F( D2F ), _G( G ), _DG( DG ), _D2G( D2G ),
                                         _maxIterations( maxIterations ), _quiet( quiet ), _eps( tolerance ),
                                         _tau( tolerance ), _n( -1 ), _fixedVariables( nullptr ) {}

  void setFixedVariables( const std::vector<int> &fixedVariables ) override {
    _fixedVariables = &fixedVariables;
  }

  void solve( const VectorType &x_0, VectorType &x_k ) const {
    this->status.Iteration = 0;
    this->status.totalTime = 0.;
    this->status.additionalTimings["TangentialStep"] = 0.;
    this->status.additionalTimings["NormalStep"] = 0.;

    // Initial value
    x_k = x_0;
    RealType trRadius = _initRadius;


    // Values at iterates
    RealType F_k;
    VectorType JacF_k, G_k, dn, dt, p_k;
    MatrixType A_k, AAT;
    VectorType lambda_k, opt_vec;
    MatrixType H_k;
    H_k.resize( x_0.size(), x_0.size());

    RealType eta_k, eps_k, rho_k;

    VectorType tmp_x_k( x_0.size()); // Temporary new iterate
    RealType tmp_F_k;
    VectorType tmp_G_k;
    RealType phi_k, tmp_phi_k; // Values of merit function
    RealType mu_k = 1., tmp_mu_k; // Penalty parameter for merit function
    RealType m_red, f_red; // predicted and actual reduction

    VectorType diagonal( x_0.size());
    MatrixType Dinv( x_0.size(), x_0.size());
    MatrixType D( x_0.size(), x_0.size());
    for ( int i = 0; i < x_0.size(); i++ ) {
      Dinv.coeffRef( i, i ) = 1;
      D.coeffRef( i, i ) = 1;
    }
    D.makeCompressed();
    Dinv.makeCompressed();
    VectorType c;
    MatrixType H, A, AATp;
    VectorType dn_p = VectorType::Zero(x_0.size());

    if ( !_quiet ) {
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Start Byrd-Omojokun SQP method with " << _maxIterations << " iterations and eps = " << _eps << "."
                << std::endl;
      std::cout << "=======================================================================" << std::endl;
//      std::cerr << "#iter -- F -- ||DF|| -- ||G||_inf -- ALF -- ||DL||_inf -- mu -- eta -- tau" << std::endl;
    }


    int iter = 0;
    while ( iter < _maxIterations ) {
      iter++;
      this->status.Iteration++;
      auto t_start = std::chrono::high_resolution_clock::now();

      // Step 0: Evaluate functionals and gradients
      _F.apply( x_k, F_k );
      _DF.apply( x_k, JacF_k );
      if ( _fixedVariables )
        applyMaskToVector( *_fixedVariables, JacF_k );

      _G.apply( x_k, G_k );
      _DG.apply( x_k, A_k );

      AAT = A_k * A_k.transpose();

      // Step 1: Prepare projections
//      if ( iter == 1 )
      _linearSolver.analyzePattern( AAT );
      _linearSolver.factorize( AAT );

      // Step 2: Compute multiplier estimates and check stopping criterion
      VectorType Ax = A_k * JacF_k;
      lambda_k = _linearSolver.solve( Ax );

      opt_vec.noalias() = JacF_k - A_k.transpose() * lambda_k;

      if ( opt_vec.template lpNorm<Eigen::Infinity>() < _eps && G_k.template lpNorm<Eigen::Infinity>() < _tau ) {
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "Acceptable solution." << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
        break;
      }

      // Step 3: Solve normal subproblem
//      SolverStatus<ConfiguratorType> normalStatus = solveNormalSubproblem( A_k, G_k, trRadius, dn );
//
//      if ( _fixedVariables )
//        applyMaskToVector( *_fixedVariables, dn );

      // Step 4: Evaluate Hessian
      evaluateHessian( x_k, lambda_k, H_k );
      if ( _fixedVariables )
        applyMaskToSymmetricMatrix( *_fixedVariables, H_k );

      // Step 5: Compute update / solve tangential subproblem
      if ( _diagonalPreconditioning == 1 ) {
        diagonal = H_k.diagonal();
        for ( int i = 0; i < x_0.size(); i++ ) {
          D.coeffRef( i, i ) = diagonal[i];
          if ( std::abs( diagonal[i] ) > eps_red )
            Dinv.coeffRef( i, i ) = 1 / diagonal[i];
          else
            Dinv.coeffRef( i, i ) = 0;
        }
      }

      SolverStatus<ConfiguratorType> tangentialStatus;
      SolverStatus<ConfiguratorType> normalStatus;

      if ( _diagonalPreconditioning == -1 ) {
        normalStatus = solveNormalSubproblem( A_k, G_k, trRadius, dn );
        tangentialStatus = solveTangentialSubproblem( H_k, A_k, JacF_k, trRadius, dn, p_k );
      }
      else {
//        VectorType dn_test = VectorType::Zero(x_0.size());
//        VectorType p_test = VectorType::Zero(x_0.size());
//        solveNormalSubproblem( A_k, G_k, trRadius, dn_p );
//        solveTangentialSubproblem( H_k, A_k, JacF_k, trRadius, dn_p, p_k );

        c.noalias() = Dinv * JacF_k;
        H = Dinv * H_k * Dinv;
        A = A_k * Dinv;
//        dn_p = D * dn;

        AATp = A * A.transpose();

//        std::cout << std::endl;
//        std::cout << std::scientific << " pre Dinv: " << std::endl
//                  << Eigen::MatrixXd( c ).format( Eigen::IOFormat( Eigen::StreamPrecision, 0, ", ", "", "", "\n", "", "" ))
//                  << std::endl;
//        std::cout << "post Dinv: " << std::endl
//                  << Eigen::MatrixXd( JacF_k ).format(Eigen::IOFormat( Eigen::StreamPrecision, 0, ", ", "", "", "\n", "", "" ))
//                  << std::endl;
//        std::cout << std::endl;

        _linearSolver.analyzePattern( AATp );
        _linearSolver.factorize( AATp );
        normalStatus = solveNormalSubproblem( A, G_k, trRadius, dn_p );
        tangentialStatus = solveTangentialSubproblem( H, A, c, trRadius, dn_p, p_k );

//        std::cout << std::endl;
//        std::cout << std::scientific << " pre Dinv: " << std::endl
//                  << Eigen::MatrixXd( p_test ).format( Eigen::IOFormat( Eigen::StreamPrecision, 0, ", ", ",", "", "", "(", ")" ))
//                  << std::endl;
//        std::cout << "post Dinv: " << std::endl
//                  << Eigen::MatrixXd( p_k ).format(Eigen::IOFormat( Eigen::StreamPrecision, 0, ", ", ",", "", "", "(", ")" ))
//                  << std::endl;
//        std::cout << std::endl;
      }

      RealType pkn = p_k.norm();

      if ( _diagonalPreconditioning != -1 ) {

        p_k = Dinv * p_k;

      }

      if ( _fixedVariables ) {
        std::cout << " fixing variables" << std::endl;
        applyMaskToVector( *_fixedVariables, p_k );
      }


      // Step 6: Evaluate model reduction and update mu (penalty parameter) if necessary
      tmp_x_k = x_k + p_k; // temporary new iterate
      _F.apply( tmp_x_k, tmp_F_k );
      _G.apply( tmp_x_k, tmp_G_k );

      std::tie( m_red, tmp_mu_k ) = evaluateModelReductionAndPenalty( p_k, G_k, JacF_k, A_k, H_k, mu_k );


      // Step 7: Compute actual improvement and quotient
      tmp_phi_k = evaluateMeritFunction( tmp_F_k, tmp_G_k, tmp_mu_k );
      phi_k = evaluateMeritFunction( F_k, G_k, tmp_mu_k );

      f_red = (phi_k - tmp_phi_k); //! \todo Reevaluate phi_k with new mu_k?

      if ((std::abs( f_red ) < eps_red && std::abs( m_red ) < eps_red) || std::abs( f_red - m_red ) < eps_red ) {
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "Cutoff active in rho-computation"
                  << std::endl;
        rho_k = 1.;
      }
      else
        rho_k = f_red / m_red; // actual over predicted reduction

      // Step X: Second-order correction
      if ( _diagonalPreconditioning == -1 &&  rho_k < _eta && dn.norm() <= _socThreshold * (p_k - dn).norm() ) {
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "Going into SOC." << std::endl;
        // Determine second-order correction
        VectorType p_soc(x_0.size());
//        std::cout << "Norm of step vs radius: " << p_k.norm() << " vs " << trRadius << std::endl;
        computeSecondOrderCorrection( H_k, A_k, JacF_k, trRadius, tmp_G_k, p_k, p_soc );

        // Potential new iterate
        VectorType x_soc =  x_k + p_soc;

        // Evaluate function and constraint
        VectorType G_soc;
        RealType F_soc;
        _F.apply( x_soc, F_soc );
        _G.apply( x_soc, G_soc );

//        std::cout << "SOC - norm of F&G: " << F_soc.norm() << " -- " << G_soc.norm() << std::endl;

        // Compute new reductions and ratio
        RealType f_red_soc =  phi_k - evaluateMeritFunction( F_soc, G_soc, tmp_mu_k );

        RealType rho_soc;
        if ((std::abs( f_red_soc ) < eps_red && std::abs( m_red ) < eps_red) || std::abs( f_red_soc - m_red ) < eps_red ) {
          std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "SOC: Cutoff active in rho-computation"
                    << std::endl;
          rho_soc = 1.;
        }
        else
          rho_soc = f_red_soc / m_red;

//        std::cout << "SOC - rho: " << rho_soc << " (" << f_red_soc << " / " << m_red << ")" << std::endl;


        if ( rho_soc > _eta ) {
          std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "SOC successful!" << std::endl;
          tmp_x_k = x_soc;
          p_k = p_soc;

          tmp_F_k = F_soc;
          tmp_G_k = G_soc;

          rho_k = rho_soc;
        }

      }

      // Step 8: Accept or decline step and update trust region radius
      if ( rho_k > _eta ) {
        x_k = tmp_x_k;
        mu_k = tmp_mu_k;
      }
      //! \todo Set flag that Hessian is not recomputed if iterate is declined

      //! \todo turn control over trust region size into parameters
      if ( rho_k < _eta ) {
        trRadius = trRadius / 4.;
      }
      else if ( rho_k > 0.75 ) {
        if ( _diagonalPreconditioning == -1 && pkn >= 0.9 * trRadius )
          trRadius = std::min( 2 * trRadius, _maxRadius );
        else if ( std::abs( pkn - trRadius ) <= eps_red )
          trRadius = std::min( 2 * trRadius, _maxRadius );
      }

      // Step 9: Additional stopping criteria
      if ( trRadius < _minRadius ) {
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "Trust region too small." << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
        break;
      }

      if ( p_k.template lpNorm<Eigen::Infinity>() < _minStepsize ) {
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << "Step size too small."
                  << " || " << pkn
                  << " ||===|| " << tangentialStatus.Residual << " || " << std::fixed << std::setw( 2 )
                  << tangentialStatus.reasonOfTermination << " || " << std::setw( 4 ) << tangentialStatus.Iteration
                  << std::endl;

        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
        break;
      }

      if ( !_quiet )
        std::cout << " -- BOSQP -- Iter " << std::setw( 3 ) << iter << ": " << std::scientific << std::setprecision( 6 )
                  << F_k
                  << " || " << G_k.template lpNorm<Eigen::Infinity>()
                  << " || " << opt_vec.template lpNorm<Eigen::Infinity>()
                  << " ||===|| " << std::setw( 13 ) << rho_k
                  << " || " << p_k.template lpNorm<Eigen::Infinity>()
                  << " || " << trRadius
                  << " || " << mu_k
                  << " ||===|| " << tangentialStatus.Residual
                  << " || " << std::fixed << std::setw( 2 ) << tangentialStatus.reasonOfTermination
                  << " || " << std::setw( 4 ) << tangentialStatus.Iteration
                  //                  << std::setprecision(2)
                  //                  << " ||===|| " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end - t_start ).count()
                  //                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_inner - t_start_inner ).count()
                  //                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_pre - t_start_pre ).count()
                  //                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_eval - t_start_eval ).count()
                  << std::endl;
    }

    if ( !_quiet ) {
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Finished Byrd-Omojokun SQP method after " << iter << " steps." << std::endl;
      std::cout << "=======================================================================" << std::endl;
      std::cout << " Iter " << std::setw( 3 ) << iter << ": " << std::scientific << std::setprecision( 6 )
                << F_k
                << " || " << G_k.template lpNorm<Eigen::Infinity>()
                << " || " << opt_vec.template lpNorm<Eigen::Infinity>()
                << " ||===|| " << std::setw( 13 ) << rho_k
                << " || " << p_k.template lpNorm<Eigen::Infinity>()
                << " || " << trRadius
                //              << std::setprecision( 2 )
                //              << " ||===|| " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end - t_start ).count()
                //              << " || " << std::setw( 6 )
                //              << std::chrono::duration<RealType, std::milli>( t_end_inner - t_start_inner ).count()
                //              << " || " << std::setw( 6 )
                //              << std::chrono::duration<RealType, std::milli>( t_end_pre - t_start_pre ).count()
                //              << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_eval - t_start_eval ).count()
                << std::endl;

      std::cout << "=======================================================================" << std::endl << std::endl;


    }
  }

protected:
  void evaluateHessian( const VectorType &x_k, const VectorType &lambda_k, MatrixType &H ) const {
    H.setZero();

    TripletListType triplets;
    _D2F.pushTriplets( x_k, triplets );
    _D2G.pushTriplets( x_k, triplets, lambda_k );

    H.setFromTriplets( triplets.begin(), triplets.end());
  }

  RealType evaluateMeritFunction( const RealType &F_k, const VectorType &G_k, const RealType &mu ) const {
    if ( _meritFunction == "L2" ) {
      // f + mu * ||c||
      return F_k + mu * G_k.template lpNorm<2>();
    }
    else
      throw std::runtime_error( "ByrdOmojokunSQP::evaluateMeritFunction(): Unknown merit function!" );
  }

  std::tuple<RealType, RealType> evaluateModelReductionAndPenalty( const VectorType &p_k, const VectorType &G_k,
                                                                   const VectorType &JacF_k, const MatrixType &A_k,
                                                                   const MatrixType &H_k, const RealType &mu ) const {
    if ( _meritFunction == "L2" ) {
      RealType quadratic_model = -JacF_k.dot( p_k ) - 0.5 * p_k.dot( H_k * p_k );
      VectorType lin_constr = G_k + A_k * p_k;
      RealType vpred = std::max( eps_red, G_k.norm() - lin_constr.norm());

      RealType new_mu = std::max( mu, quadratic_model / ((_penaltyFactor - 1.) * vpred));

      RealType mred = new_mu * vpred + quadratic_model;


      return std::make_tuple( mred, new_mu );
    }
    else
      throw std::runtime_error( "ByrdOmojokunSQP::evaluateModelReduction(): Unknown merit function!" );
  }

  SolverStatus<ConfiguratorType> solveNormalSubproblem( const MatrixType &A_k, const VectorType &G_k,
                                                        const RealType &trRadius, VectorType &v ) const {
    if ( _normalMethod == "ModifiedDogleg" ) {
      RealType reducedRadius = 0.8 * trRadius;

      ModifiedDoglegSolver<ConfiguratorType, LinearSolverType> modogSolver( A_k, G_k, reducedRadius, _linearSolver );
      modogSolver.setParameters( _normalIntParameters );
      modogSolver.setParameters( _normalRealParameters );
      modogSolver.setParameters( _normalStringParameters );
      modogSolver.solve( v, v );
      return modogSolver.Status();
    }

    else
      throw std::runtime_error( "ByrdOmojokunSQP::solveNormalSubproblem(): Unknown normal subproblem method!" );
  }

  SolverStatus<ConfiguratorType> solveTangentialSubproblem( const MatrixType &H_k, const MatrixType &A_k,
                                                            const VectorType &JacF_k, const RealType &radius,
                                                            const VectorType &n, VectorType &p ) const {
    if ( _tangentialMethod == "ProjectedCG" ) {
      RealType reducedRadius;
      if (_diagonalPreconditioning == -1)
        reducedRadius = 0.9 * radius;
      else
        reducedRadius = radius;

      VectorType b( A_k * n );

      ProjectedCGSolver<ConfiguratorType, LinearSolverType> prCGSolver( H_k, JacF_k, A_k, b, reducedRadius, _linearSolver,
                                                                        1e-8, n.size(), true );

      prCGSolver.setParameters( _tangentialIntParameters );
      prCGSolver.setParameters( _tangentialRealParameters );
      prCGSolver.setParameters( _tangentialStringParameters );

      prCGSolver.solve( n, p );

      return prCGSolver.Status();
    }
    else
      throw std::runtime_error( "ByrdOmojokunSQP::solveTangentialSubproblem(): Unknown tangential subproblem method!" );
  }

  void computeSecondOrderCorrection( const MatrixType &H_k, const MatrixType &A_k,
                                                               const VectorType &JacF_k, const RealType &radius,
                                                               const VectorType &G_p,
                                                               const VectorType &p, VectorType &p_soc ) const {
//      VectorType lSol = _linearSolver.solve(G_p);
//      p_soc = p + A_k.transpose() * lSol;


      VectorType b( A_k * p - G_p );

      ProjectedCGSolver<ConfiguratorType, LinearSolverType> prCGSolver( H_k, JacF_k, A_k, b, radius, _linearSolver,
                                                                        1e-8, p.size(), true );

      prCGSolver.setParameters( _tangentialIntParameters );
      prCGSolver.setParameters( _tangentialRealParameters );
      prCGSolver.setParameters( _tangentialStringParameters );

      prCGSolver.solve( p, p_soc );

  }


public:
  void setParameter( const std::string &name, std::string value ) override {
    if ( name == "tangential_method" )
      _tangentialMethod = value; //! \todo Add safeguard against undefined methods
    else if ( name == "normal_method" )
      _normalMethod = value; //! \todo Add safeguard against undefined methods
    else if ( name.rfind( "normal__", 0 ) == 0 )
      _normalStringParameters[name.substr( 8, std::string::npos )] = value;
    else if ( name.rfind( "tangential__", 0 ) == 0 )
      _tangentialStringParameters[name.substr( 12, std::string::npos )] = value;
    else
      throw std::logic_error( "ByrdOmojokunSQP::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "constraint_tolerance" )
      _tau = value;
    else if ( name == "optimality_tolerance" )
      _eps = value;
    else if ( name == "maximal_radius" )
      _maxRadius = value;
    else if ( name == "initial_radius" )
      _initRadius = value;
    else if ( name == "minimal_radius" )
      _minRadius = value;
    else if ( name == "minimal_stepsize" )
      _minStepsize = value;
    else if ( name == "soc_threshold" )
      _socThreshold = value;
    else if ( name.rfind( "normal__", 0 ) == 0 )
      _normalRealParameters[name.substr( 8, std::string::npos )] = value;
    else if ( name.rfind( "tangential__", 0 ) == 0 )
      _tangentialRealParameters[name.substr( 12, std::string::npos )] = value;
    else
      throw std::runtime_error( "ByrdOmojokunSQP::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else if ( name == "diagonal_preconditioning" )
      _diagonalPreconditioning = value;
    else if ( name.rfind( "normal__", 0 ) == 0 )
      _normalIntParameters[name.substr( 8, std::string::npos )] = value;
    else if ( name.rfind( "tangential__", 0 ) == 0 )
      _tangentialIntParameters[name.substr( 12, std::string::npos )] = value;
    else
      throw std::runtime_error( "ByrdOmojokunSQP::setParameter(): Unknown parameter '" + name + "'." );
  }

};

#endif //OPTIMIZATION_BYRDOMOJOKUNSQP_H
