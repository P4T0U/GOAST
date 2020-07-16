// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __CONFIGURATORS_H
#define __CONFIGURATORS_H

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC system_header
#endif


//== INCLUDES =================================================================
#include "SmallVecMat.h"
#include "Auxiliary.h"


//==========================================================================================================
// DEFAULT CONFIGURATOR
//==========================================================================================================
struct DefaultConfigurator {
public:
  // Underlying floating-point type
  typedef double RealType;

  // Dense n-dimensional linear algebra
  typedef Eigen::VectorXd VectorType;
  typedef Eigen::MatrixXd FullMatrixType;

  // Sparse linear algebra
  // second template argument of SparseMatrix: union of bit flags controlling the storage scheme. Currently the only possibility is ColMajor or RowMajor. The default is 0 which means column-major.
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> SparseMatrixType;
  typedef Eigen::Triplet<double> TripletType;
  typedef GenericTensor<SparseMatrixType> TensorType;

  // 3-dimensional linear algebra
  typedef SmallVec3<RealType> VecType;
  typedef SmallMat33<RealType> MatType;

  // Rotations and Orientations
  typedef Eigen::Matrix<RealType, 3, 3> RotationType;
  typedef Eigen::Matrix<RealType, 3, 3> FrameType;
  typedef Eigen::Quaternion<RealType> QuaternionType;
  typedef Eigen::AngleAxis<RealType> AngleAxisType;
};

#endif