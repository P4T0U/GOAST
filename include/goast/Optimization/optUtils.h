// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Different utility functions for optimization algorithms
 * \author Sassen
 */

#ifndef OPTIMIZATION_OPTUTILS_H
#define OPTIMIZATION_OPTUTILS_H

#include <cholmod.h>

#include <tuple>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * \brief Compute the segment of a line (segment) u+t*v intersecting with a coordinate box
 * \author Sassen
 * \param u base point of the line
 * \param v direction vector of the line
 * \param lowerBounds lower bounds of the box for the different coordinates
 * \param upperBounds upper bounds of the box for the different coordinates
 * \param min_t Minimal value of t for the segment of the line to be considered
 * \param max_t Maximal value of t for the segment of the line to be considered
 * \return Tuple with the starting and ending t of the intersecting segment and a bool indicating if there is an intersection
 */
template<typename RealType, typename VectorType>
std::tuple<RealType, RealType, bool> lineBoxIntersection( const VectorType &u, const VectorType &v,
                                                          const VectorType &lowerBounds, const VectorType &upperBounds,
                                                          const RealType min_t = 0., const RealType max_t = 1. ) {
  if ( v.norm() == 0 )
    return std::make_tuple( 0., 0., false );

  int last_index = 0;
  RealType r_1 = -std::numeric_limits<RealType>::infinity();
  RealType r_2 = std::numeric_limits<RealType>::infinity();

  for ( int i = 0; i < v.size(); i++ ) {
    // Check if the boundaries are not satisfied for some coordinates for which v is zero,
    // then there is no intersection
    if ( v[i] == 0 ) {
      if ( u[i] < lowerBounds[last_index] || u[i] > upperBounds[last_index] )
        return std::make_tuple( 0., 0., false );
    }
    else {
      // For non-zero entries
      RealType z_value = u.coeff( i );
      RealType r_lb = (lowerBounds[i] - z_value) / v[i];
      RealType r_ub = (upperBounds[i] - z_value) / v[i];
      r_1 = std::max( r_1, std::min( r_lb, r_ub ));
      r_2 = std::min( r_2, std::max( r_lb, r_ub ));
    }
  }

  if ( r_1 > r_2 )
    return std::make_tuple( 0., 0., false );

  if ( r_2 < min_t || r_1 > max_t )
    return std::make_tuple( 0., 0., false );
  else
    return std::make_tuple( std::max( min_t, r_1 ), std::min( max_t, r_2 ), true );
}

/**
 * \brief Compute the segment of a line (segment) u+t*v intersecting with a ball
 * \author Sassen
 * \param u base point of the line
 * \param v direction vector of the line
 * \param radius radius of the ball
 * \param min_t Minimal value of t for the segment of the line to be considered
 * \param max_t Maximal value of t for the segment of the line to be considered
 * \return Tuple with the starting and ending t of the intersecting segment and a bool indicating if there is an intersection
 */
template<typename RealType, typename VectorType>
std::tuple<RealType, RealType, bool> lineBallIntersection( const VectorType &u,
                                                           const VectorType &v,
                                                           const RealType radius,
                                                           const RealType min_t = -std::numeric_limits<RealType>::infinity(),
                                                           const RealType max_t = std::numeric_limits<RealType>::infinity()) {
  if ( v.norm() == 0 )
    return std::make_tuple( 0., 0., false );

  if ( radius == std::numeric_limits<RealType>::infinity())
    return std::make_tuple( min_t, max_t, true );

  // coefficients of quadratic equation
  RealType a = v.squaredNorm();
  RealType b = 2 * u.dot( v );
  RealType c = u.squaredNorm() - radius * radius;

  // abc-Formula
  RealType discriminant = b * b - 4 * a * c;

  if ( discriminant < 0 )
    return std::make_tuple( 0., 0., false );

  // For stability reasons we first determine the root where b and the discriminant have the same sign
  // and then derive the by the relation r_2 = (-c/a)/r_1
  RealType aux = b + std::copysign( std::sqrt( discriminant ), b );
  RealType r_1 = -aux / (2 * a);
  RealType r_2 = -2 * c / aux;

  RealType r_a, r_b;
  std::tie( r_a, r_b ) = std::minmax( r_1, r_2 );

  if ( r_b < min_t || r_a > max_t )
    return std::make_tuple( 0., 0., false );
  else
    return std::make_tuple( std::max( min_t, r_a ), std::min( max_t, r_b ), true );
}

/**
 * \brief Compute the segment of a line u+t*v intersecting with a sphere
 * \author Sassen
 * \param u base point of the line
 * \param v direction vector of the line
 * \param radius radius of the sphere
 * \return Tuple with the starting and ending t of the intersecting segment and a bool indicating if there is an intersection
 */
template<typename RealType, typename VectorType>
std::tuple<RealType, RealType, bool> lineSphereIntersection( const VectorType &u,
                                                             const VectorType &v,
                                                             const RealType radius ) {
  if ( v.norm() == 0 )
    return std::make_tuple( 0., 0., false );

  if ( radius == std::numeric_limits<RealType>::infinity())
    throw std::domain_error( "Line-Sphere-Intersection not defined for sphere with infinite radius." );

  // coefficients of quadratic equation
  RealType a = v.squaredNorm();
  RealType b = 2 * u.dot( v );
  RealType c = u.squaredNorm() - radius * radius;

  // abc-Formula
  RealType discriminant = b * b - 4 * a * c;

  if ( discriminant < 0 )
    return std::make_tuple( 0., 0., false );


  // For stability reasons we first determine the root where b and the discriminant have the same sign
  // and then derive the by the relation r_2 = (-c/a)/r_1
  RealType aux = b + std::copysign( std::sqrt( discriminant ), b );
  RealType r_1 = -aux / (2 * a);
  RealType r_2 = -2 * c / aux;

  RealType r_a, r_b;
  std::tie( r_a, r_b ) = std::minmax( r_1, r_2 );

  return std::make_tuple( r_a, r_b, true );
}

/**
 * \brief Compute the segment of a line (segment) u+t*v intersecting with a coordinate box and a ball
 * \author Sassen
 * \param u base point of the line
 * \param v direction vector of the line
 * \param radius radius of the ball
 * \param lowerBounds lower bounds of the box for the different coordinates
 * \param upperBounds upper bounds of the box for the different coordinates
 * \param min_t Minimal value of t for the segment of the line to be considered
 * \param max_t Maximal value of t for the segment of the line to be considered
 * \return Tuple with the starting and ending t of the intersecting segment and a bool indicating if there is an intersection
 */
template<typename RealType, typename VectorType>
std::tuple<RealType, RealType, bool> lineBoxBallIntersections( const VectorType &u, const VectorType &v,
                                                               RealType radius,
                                                               const VectorType &lowerBounds,
                                                               const VectorType &upperBounds,
                                                               const RealType min_t = 0.,
                                                               const RealType max_t = 1. ) {
  RealType sph_r_1, sph_r_2, box_r_1, box_r_2;
  bool sph_intersect, box_intersect;

  std::tie( sph_r_1, sph_r_2, sph_intersect ) = lineBallIntersection( u, v, radius, min_t, max_t );
  std::tie( box_r_1, box_r_2, box_intersect ) = lineBoxIntersection( u, v, lowerBounds, upperBounds, min_t, max_t );

  if ( sph_intersect && box_intersect ) {
    RealType r_1 = std::max( sph_r_1, box_r_1 );
    RealType r_2 = std::min( sph_r_2, box_r_2 );
    if ( r_1 <= r_2 )
      return std::make_tuple( r_1, r_2, true );
  }

  return std::make_tuple( 0., 0., false );
}

/**
 * \brief Check if a point is inside coordinate box
 * \param x point
 * \param lowerBounds lower bounds of the box for the different coordinates
 * \param upperBounds upper bounds of the box for the different coordinates
 * \return whether the point is inside the coordinate box or not
 */
template<typename VectorType>
bool insideBox( const VectorType &x, const VectorType &lowerBounds, const VectorType &upperBounds ) {
  int last_index = 0;
  for ( int i = 0; i < x.size(); i++ ) {
    if ( x[i] < lowerBounds[i] || x[i] > upperBounds[i] )
      return false;
  }

  return true;
}


template<typename Scalar, int Flags, typename StorageIndex>
Eigen::Map<Eigen::SparseMatrix<Scalar, Flags, StorageIndex>> viewAsEigen( cholmod_factor &cf ) {
  return Eigen::Map<Eigen::SparseMatrix<Scalar, Flags, StorageIndex>>
          ( cf.n, cf.n, static_cast<StorageIndex *>(cf.p)[cf.n],
            static_cast<StorageIndex *>(cf.p), static_cast<StorageIndex *>(cf.i), static_cast<Scalar *>(cf.x));
}

#endif //OPTIMIZATION_OPTUTILS_H
