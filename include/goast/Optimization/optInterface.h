// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for defining the abstract interfaces to optimization methods
 */

#ifndef OPTIMIZATION_INTERFACE_H
#define OPTIMIZATION_INTERFACE_H

#include <goast/Core/Auxiliary.h>

#include "Objectives.h"
#include "Constraints.h"
#include "optUtils.h"

template<typename ConfiguratorType>
struct SolverStatus {
  int Iteration;
  typename ConfiguratorType::RealType Residual;
  typename ConfiguratorType::RealType totalTime;
  std::map<std::string, typename ConfiguratorType::RealType> additionalTimings;
  std::map<std::string, int> additionalIterations;
  int reasonOfTermination;
};

template<typename ConfiguratorType>
class OptimizationBase {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

public:
  /**
   * \brief Status of the solver
   */
  mutable SolverStatus<ConfiguratorType> status;

  OptimizationBase() = default;

  virtual ~OptimizationBase() = default;

  const SolverStatus<ConfiguratorType> &Status() const {
    return status;
  }


  /**** Mandatory methods for solving the problem ****/

  /**
   * \brief Solve the optimization problem with given starting point
   * \param[in] startingPoint given starting point
   * \param[out] solution the computed solution
   */
  virtual void solve( const VectorType &startingPoint, VectorType &solution ) const = 0;


  /**** Mandatory methods for setting parameters ****/

  /**
   * \brief Set a real-valued parameter
   * \param[in] name Name of the parameter
   * \param[in] value Value of the parameter
   */
  virtual void setParameter( const std::string &name, RealType value ) = 0;

  /**
   * \brief Set a integer-valued parameter
   * \param[in] name Name of the parameter
   * \param[in] value Value of the parameter
   */
  virtual void setParameter( const std::string &name, int value ) = 0;

  /**
   * \brief Set a string-valued parameter
   * \param[in] name Name of the parameter
   * \param[in] value Value of the parameter
   */
  virtual void setParameter( const std::string &name, std::string value ) = 0;

  /**** Optional methods for defining the problem ****/

  /**
   * \brief Set which variables should remain fixed during the optimization
   * \param[in] fixedVariables Vector containing indices of variables to fix
   * \note This is an optional method, not every method has to implement it.
   */
  virtual void setFixedVariables( const std::vector<int> &fixedVariables ) {
    throw std::logic_error( "OptimizationBase::setFixedVariables(): Unimplemented function! The chosen method does "
                            "not seem to provide the possibility to fix variables." );
  }

  /**
   * \brief Set box constraints i.e. lower and upper bounds on the individual variables
   * \param[in] lowerBounds Vector containing lower bounds on the variables
   * \param[in] upperBounds Vector containing upper bounds on the variables
   */
  virtual void setVariableBounds( const VectorType &lowerBounds, const VectorType &upperBounds ) {
    throw std::logic_error( "OptimizationBase::setVariableBounds(): Unimplemented function! The chosen method does "
                            "not seem to allow box constraints." );
  }

  /**** Helper functions ****/
  /**
   * \brief Set int-valued parameters in batch
   * \param Parameters map containing parameters to set
   */
  void setParameters( const std::map<std::string, int> &Parameters ) {
    for (auto const& p : Parameters )
      setParameter(p.first, p.second);
  }

  /**
   * \brief Set string-valued parameters in batch
   * \param Parameters map containing parameters to set
   */
  void setParameters( const std::map<std::string, std::string> &Parameters ) {
    for (auto const& p : Parameters )
      setParameter(p.first, p.second);
  }

  /**
   * \brief Set real-valued parameters in batch
   * \param Parameters map containing parameters to set
   */
  void setParameters( const std::map<std::string, RealType> &Parameters ) {
    for (auto const& p : Parameters )
      setParameter(p.first, p.second);
  }
};

#endif //OPTIMIZATION_INTERFACE_H
