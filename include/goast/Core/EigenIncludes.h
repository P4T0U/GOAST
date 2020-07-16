// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __EIGENINCLUDES_H
#define __EIGENINCLUDES_H

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC system_header
#endif


//== INCLUDES =================================================================

#ifdef GOAST_WITH_MKL
#define EIGEN_USE_MKL_ALL
#endif

// EIGEN
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>

// Custom adaptions, e.g. reductions, for Eigen and OpenMP
#ifdef _OPENMP
#pragma omp declare reduction (+: Eigen::VectorXd: omp_out=omp_out+omp_in) initializer(omp_priv=Eigen::VectorXd::Zero(omp_orig.size()))
#endif

#endif