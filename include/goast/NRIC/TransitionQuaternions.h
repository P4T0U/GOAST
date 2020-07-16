// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Implementation of intrinsic transition rotations as quaternions
 * \author Sassen
 *
 * Based on Yue Wang, B. Liu, and Y. Tong. Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Comput. Graph. Forum, 31:2277–2287, 2012.
 *
 * \note This directly relies on the Axis-Angle and Quaternion implementations from Eigen
 */

#ifndef NRIC_TRANSITIONQUATERNIONS_H
#define NRIC_TRANSITIONQUATERNIONS_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/Topology.h>

#include "TrigonometryOperators.h"

#include <algorithm>
#include <array>

#include <Eigen/Geometry>

/**
 * \brief Compute local transition quaternion from edge lengths and dihedral angle
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This computes the rotation \f$Q(\theta,a,b,c) = R_x(\theta) R_z(\gamma(a,b,c))\f$, where \gamma(a,b,c) is the
 * interior angle determined via law of cosines
 */
template<typename ConfiguratorType>
class LocalTransitionQuaternionOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::AngleAxisType AngleAxisType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef typename ConfiguratorType::QuaternionType RotationType;
  const Eigen::Matrix<RealType, 3, 1> UnitZ = Eigen::Matrix<RealType, 3, 1>::UnitZ();
  const Eigen::Matrix<RealType, 3, 1> UnitX = Eigen::Matrix<RealType, 3, 1>::UnitX();

public:
  explicit LocalTransitionQuaternionOp() = default;

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    const RealType &theta = Arg[0];
    const RealType &a = Arg[1];
    const RealType &b = Arg[2];
    const RealType &c = Arg[3];

    RealType gamma = acos((pow( a, 2 ) + pow( b, 2 ) - pow( c, 2 )) / (2 * a * b));
    QuaternionType q( AngleAxisType( theta, UnitX ) * AngleAxisType( gamma, UnitZ ));

    Dest.resize( 4 );
    Dest.setZero();

    Dest[0] = q.w();
    Dest[1] = q.x();
    Dest[2] = q.y();
    Dest[3] = q.z();
  }

  void apply( const RealType &theta, const RealType &a, const RealType &b, const RealType &c,
              RotationType &Dest ) const {
    RealType gamma = acos((pow( a, 2 ) + pow( b, 2 ) - pow( c, 2 )) / (2 * a * b));
    Dest = AngleAxisType( theta, UnitX ) * AngleAxisType( gamma, UnitZ );
  }

//  RotationType operator()( const RealType &theta, const RealType &a, const RealType &b, const RealType &c ) {
//    RealType gamma = acos((pow( a, 2 ) + pow( b, 2 ) - pow( c, 2 )) / (2 * a * b));
//    return RotationType( AngleAxisType( theta, UnitX ) * AngleAxisType( gamma, UnitZ ));
//  }

};

/**
 * \brief Compute derivative of local transition quaternion w.r.t. edge lengths and dihedral angle
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see LocalTransitionQuaternionOp
 */
template<typename ConfiguratorType>
class LocalTransitionQuaternionGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::AngleAxisType AngleAxisType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef typename ConfiguratorType::QuaternionType RotationType;
  const Eigen::Matrix<RealType, 3, 1> UnitZ = Eigen::Matrix<RealType, 3, 1>::UnitZ();
  const Eigen::Matrix<RealType, 3, 1> UnitX = Eigen::Matrix<RealType, 3, 1>::UnitX();


public:
  explicit LocalTransitionQuaternionGradient() = default;

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    std::array<QuaternionType, 4> qdiff;
    apply( Arg[0], Arg[1], Arg[2], Arg[3], qdiff );

    Dest.resize( 4, 4 );
    Dest.setZero();

    for ( int i = 0; i < 4; i++ ) {
      Dest.coeffRef( 0, i ) = qdiff[i].w();
      Dest.coeffRef( 1, i ) = qdiff[i].x();
      Dest.coeffRef( 2, i ) = qdiff[i].y();
      Dest.coeffRef( 3, i ) = qdiff[i].z();
    }
  }

  void apply( const RealType &theta, const RealType &a, const RealType &b, const RealType &c,
              std::array<QuaternionType, 4> &Dest ) const {
    //! \todo Use numerically stable schemes for this
    const RealType gamma = std::acos((a * a + b * b - c * c) / (2 * a * b));
    const RealType Area = std::sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))) / 4.;

    // Evaluate trigonometric functions
    const RealType sinThetaHalf = std::sin( theta / 2 );
    const RealType cosThetaHalf = std::cos( theta / 2 );
    const RealType sinGammaHalf = std::sin( gamma / 2 );
    const RealType cosGammaHalf = std::cos( gamma / 2 );

    // Trigonometric products
    const RealType sinTcosG = sinThetaHalf * cosGammaHalf;
    const RealType cosTcosG = cosThetaHalf * cosGammaHalf;
    const RealType cosTsinG = cosThetaHalf * sinGammaHalf;
    const RealType sinTsinG = sinThetaHalf * sinGammaHalf;

    // S polynomials
    const RealType S_1 = -a * a + b * b + c * c;
    const RealType S_2 = a * a - b * b + c * c;

    QuaternionType q;
    // Dtheta
    q = QuaternionType( -sinTcosG, cosTcosG, -cosTsinG, -sinTsinG );
    q.coeffs() *= (1 / 2.);

    Dest[0] = QuaternionType( -sinTcosG, cosTcosG, -cosTsinG, -sinTsinG );
    Dest[0].coeffs() *= (1 / 2.);

    // Da
    q = QuaternionType( cosTsinG, sinTsinG, sinTcosG, -cosTcosG );
    q.coeffs() *= (S_2 / (8. * a * Area));

    Dest[1] = q;

    // Db
    q = QuaternionType( cosTsinG, sinTsinG, sinTcosG, -cosTcosG );
    q.coeffs() *= (S_1 / (8. * b * Area));
    Dest[2] = q;

    // Dc
    q = QuaternionType( -cosTsinG, -sinTsinG, -sinTcosG, cosTcosG );
    q.coeffs() *= (c / (4 * Area));
    Dest[3] = q;
  }

  std::array<QuaternionType, 4> operator()( const RealType &theta, const RealType &a, const RealType &b,
                                            const RealType &c ) {
    std::array<QuaternionType, 4> Dest;
    apply( theta, a, b, c, Dest );
    return Dest;
  }

};

/**
 * \brief Compute Hessian of transition quaternion w.r.t. edge lengths and dihedral angle
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see LocalTransitionQuaternionOp
 */
template<typename ConfiguratorType>
class LocalTransitionQuaternionHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef typename ConfiguratorType::QuaternionType RotationType;
  typedef typename ConfiguratorType::TensorType TensorType;


public:
  explicit LocalTransitionQuaternionHessian() = default;

  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    std::array<std::array<QuaternionType, 4>, 4> qdiff;
    apply( Arg[0], Arg[1], Arg[2], Arg[3], qdiff );

    Dest.resize( 4, 4, 4 );
    Dest.setZero();

    for ( int i = 0; i < 4; i++ ) {
      for ( int j = 0; j < 4; j++ ) {
        Dest[0].coeffRef( i, j ) = qdiff[i][j].w();
        Dest[1].coeffRef( i, j ) = qdiff[i][j].x();
        Dest[2].coeffRef( i, j ) = qdiff[i][j].y();
        Dest[3].coeffRef( i, j ) = qdiff[i][j].z();
      }
    }
  }

  void apply( const RealType &theta, const RealType &a, const RealType &b, const RealType &c,
              std::array<std::array<QuaternionType, 4>, 4> &Dest ) const {
    for ( int i = 0; i < 4; i++ )
      for ( int j = 0; j < 4; j++ )
        Dest[i][j] = QuaternionType( 0, 0, 0, 0 );

    //! \todo Use numerically stable schemes for this
    const RealType gamma = std::acos((a * a + b * b - c * c) / (2 * a * b));
    const RealType Area = std::sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))) / 4.;

    // Squared lengths
    const RealType aSqr = a * a;
    const RealType bSqr = b * b;
    const RealType cSqr = c * c;

    const RealType areaSqr16 = 16 * Area * Area;
    const RealType areaCub64 = 64 * Area * Area * Area;

    // Evaluate trigonometric functions
    const RealType sinThetaHalf = std::sin( theta / 2 );
    const RealType cosThetaHalf = std::cos( theta / 2 );
    const RealType sinGammaHalf = std::sin( gamma / 2 );
    const RealType cosGammaHalf = std::cos( gamma / 2 );

    // Trigonometric products
    const RealType sinTcosG = sinThetaHalf * cosGammaHalf;
    const RealType cosTcosG = cosThetaHalf * cosGammaHalf;
    const RealType cosTsinG = cosThetaHalf * sinGammaHalf;
    const RealType sinTsinG = sinThetaHalf * sinGammaHalf;

    // S polynomials
    const RealType S_1 = -a * a + b * b + c * c;
    const RealType S_2 = a * a - b * b + c * c;
    const RealType S_3 = aSqr + bSqr - cSqr;

    QuaternionType q;
    // Dtheta Dtheta
    q = QuaternionType( -cosTcosG, -sinTcosG, sinTsinG, -cosTsinG );
    q.coeffs() *= (1 / 4.);

    Dest[0][0] = q;

    // Da Dtheta
    q = QuaternionType( -sinTsinG, cosTsinG, cosTcosG, sinTcosG );
    q.coeffs() *= (S_2 / (16. * a * Area));

    Dest[0][1] = q;
    Dest[1][0] = q;

    // Db Dtheta
    q = QuaternionType( -sinTsinG, cosTsinG, cosTcosG, sinTcosG );
    q.coeffs() *= (S_1 / (16. * b * Area));
    Dest[2][0] = q;
    Dest[0][2] = q;

    // Dc Dtheta
    q = QuaternionType( sinTsinG, -cosTsinG, -cosTcosG, -sinTcosG );
    q.coeffs() *= (c / (8 * Area));
    Dest[3][0] = q;
    Dest[0][3] = q;

    // Da Da
    RealType term1 = 2. * S_3 * S_2 * S_2 * sinGammaHalf
                     - 4. * Area * S_2 * S_2 * cosGammaHalf
                     + 4. * areaSqr16 * (bSqr - cSqr) * sinGammaHalf;
    RealType term2 = term1;
    RealType term3 = 2. * S_3 * S_2 * S_2 * cosGammaHalf
                     + 4. * Area * S_2 * S_2 * sinGammaHalf
                     + 4. * areaSqr16 * (bSqr - cSqr) * cosGammaHalf;
    RealType term4 = -term3;
    q = QuaternionType( cosThetaHalf * term1, sinThetaHalf * term2, sinThetaHalf * term3, cosThetaHalf * term4 );
    q.coeffs() *= 1 / (4. * aSqr * areaCub64);

    Dest[1][1] = q;

    // Db Da
    term1 = -4. * aSqr * bSqr * cSqr * sinGammaHalf - Area * S_1 * S_2 * cosGammaHalf;
    term2 = term1;
    term3 = 4. * aSqr * bSqr * cSqr * cosGammaHalf - Area * S_1 * S_2 * sinGammaHalf;
    term4 = term3;
    q = QuaternionType( cosThetaHalf * term1, sinThetaHalf * term2, -sinThetaHalf * term3, cosThetaHalf * term4 );
    q.coeffs() *= 1 / (a * b * areaCub64);

    Dest[1][2] = q;
    Dest[2][1] = q;

    // Dc Da
    term1 = aSqr * S_1 * sinGammaHalf + Area * S_2 * cosGammaHalf;
    term2 = term1;
    term3 = -aSqr * S_1 * cosGammaHalf + Area * S_2 * sinGammaHalf;
    term4 = term3;
    q = QuaternionType( cosThetaHalf * term1, sinThetaHalf * term2, -sinThetaHalf * term3, cosThetaHalf * term4 );
    q.coeffs() *= 2 * c / (a * areaCub64);

    Dest[1][3] = q;
    Dest[3][1] = q;

    // Db Db
    term1 = 2. * S_3 * S_1 * S_1 * sinGammaHalf
            - 4. * Area * S_1 * S_1 * cosGammaHalf
            + 4. * areaSqr16 * (aSqr - cSqr) * sinGammaHalf;
    term2 = term1;
    term3 = 2. * S_3 * S_1 * S_1 * cosGammaHalf
            + 4. * Area * S_1 * S_1 * sinGammaHalf
            + 4. * areaSqr16 * (aSqr - cSqr) * cosGammaHalf;
    term4 = -term3;
    q = QuaternionType( cosThetaHalf * term1, sinThetaHalf * term2, sinThetaHalf * term3, cosThetaHalf * term4 );
    q.coeffs() *= 1 / (4. * bSqr * areaCub64);

    Dest[2][2] = q;

    // Dc Db
    term1 = bSqr * S_2 * sinGammaHalf + Area * S_1 * cosGammaHalf;
    term2 = term1;
    term3 = -bSqr * S_2 * cosGammaHalf + Area * S_1 * sinGammaHalf;
    term4 = term3;
    q = QuaternionType( cosThetaHalf * term1, sinThetaHalf * term2, -sinThetaHalf * term3, cosThetaHalf * term4 );
    q.coeffs() *= 2 * c / (b * areaCub64);

    Dest[2][3] = q;
    Dest[3][2] = q;

    // Dc Dc
    term1 = S_1 * S_2 * sinGammaHalf + 4 * cSqr * Area * cosGammaHalf;
    term2 = term1;
    term3 = -S_1 * S_2 * cosGammaHalf + 4 * cSqr * Area * sinGammaHalf;
    term4 = term3;
    q = QuaternionType( -cosThetaHalf * term1, -sinThetaHalf * term2, sinThetaHalf * term3, -cosThetaHalf * term4 );
    q.coeffs() *= 1 / areaCub64;

    Dest[3][3] = q;
    Dest[3][3] = q;
  }


  std::array<std::array<QuaternionType, 4>, 4> operator()( const RealType &theta, const RealType &a, const RealType &b,
                                                           const RealType &c ) {
    std::array<std::array<QuaternionType, 4>, 4> Dest;
    apply( theta, a, b, c, Dest );
    return Dest;
  }

};

#endif //NRIC_TRANSITIONQUATERNIONS_H
